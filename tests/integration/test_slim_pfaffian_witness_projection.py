"""TDD-red: assert slim Pfaffian witness projects onto conduction bands,
not band-1 sites. Expected to FAIL until A.2 lands.

Pins the slim-Pfaffian row-index bug at
`src/physics/topological_analysis.f90:1722, :1827`: the witness uses
`idx = [7, 8, n_sp+7, n_sp+8]` which assumes row `r` corresponds to band
`r`, but the wire H0 uses band-major layout (`hamiltonian_wire.f90:1208-1209`,
`row = (band-1)*N + site`). For any wire with N>1 this reads rows 7,8 within
band 1 (HH-up) instead of bands 7-8 (CB). The post-fix behavior (multi-site
scan, band-major) yields Pf != 0 at the topological-phase grid points.

For TDD-red, `call_witness_slim` is a stub that simulates the buggy
implementation by returning 0 (the band-1 Pf). The test asserts at least
one grid point returns Pf != 0; once A.2 lands and the stub is replaced by
a real `wire_pfaffian_witness_sweep` invocation, this test passes.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
# Canonical golden reference (5 B x 2 mu = 10 rows of (B, mu, z2, gap))
# produced by `build/src/topologicalAnalysis tests/regression/configs/
# wire_inas_gaas_bdg_topological_2d.toml`. Note: the brief text referenced
# `output/wire_inas_gaas_bdg_topological_2d_phase.dat`, but the executable
# actually emits `output/z2_phase_diagram.dat` and the golden lives under
# tests/regression/data/. We consume the golden for determinism (in source
# control, immune to executable drift on unrelated changes).
PHASE_FILE = (
    REPO / "tests" / "regression" / "data"
    / "wire_inas_gaas_bdg_topological_2d_phase.dat"
)


def load_phase_grid(path):
    """Read the 5x2 phase grid from the regression golden.

    Format per `topological_analysis.f90` emitter: `# B(T) mu(eV) z2 gap(eV)`.
    """
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            rows.append((float(parts[0]), float(parts[1])))
    return rows


def call_witness_slim(B, mu, N=3):
    """Call wire_pfaffian_witness_sweep at one (B, mu) point. Returns Pf sign.

    TDD-red stub: returns 0 to simulate the buggy band-1 read (rows 7,8 of
    band 1 -- which has no pairing structure at this site set -- gives
    Pf=0). After A.2 lands, this is replaced by a real Fortran invocation
    that calls wire_pfaffian_witness_sweep and returns Pf_sign.
    """
    return 0


def main():
    if not PHASE_FILE.exists():
        print(f"FAIL: {PHASE_FILE} not found")
        sys.exit(1)
    grid = load_phase_grid(PHASE_FILE)
    non_zero_count = 0
    band1_count = 0
    conduction_count = 0
    for B, mu in grid:
        pf = call_witness_slim(B, mu, N=3)
        if pf != 0:
            non_zero_count += 1
        band1_count += 1 if pf == 0 else 0
        conduction_count += 1 if pf != 0 else 0
    if non_zero_count == 0:
        print(
            f"FAIL: slim Pfaffian witness returned 0 across all {len(grid)} "
            "grid points"
        )
        print(
            f"  ({band1_count} band-1 reads, {conduction_count} "
            "conduction reads)"
        )
        print(
            "  This indicates the witness reads rows 7,8 (band-1) "
            "instead of bands 7-8 (CB)"
        )
        sys.exit(1)
    print(
        f"PASS: slim Pfaffian witness non-zero at {non_zero_count}/"
        f"{len(grid)} points"
    )


if __name__ == "__main__":
    main()